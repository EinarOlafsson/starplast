"""Strategies 33 and 34: one integrated neighbour space, and the edges the networks are missing.

Two additions to the catalogue, both built on `graphspace`, both in the eighth family ("Combine
strategies") because that is their mechanism: every other network strategy reads ONE layer, and these
two put every permitted layer and the measurement table into a single graph.

    33  Read a gene's neighbourhood in the integrated space
    34  Train on the networks and rank the edges they are missing

They share a self-test, and it is not the one strategy 16 uses. Strategy 16 scores held-out crosslinks
against RANDOM non-pairs and reports AUROC 0.97; most of that number is the discovery that
well-connected genes are well connected, since a random pair of genes is a pair of obscure genes.
Here the hidden edges are scored against DEGREE-MATCHED non-pairs -- a non-pair with a gene of similar
connectivity at both ends -- and the null re-run for the bar is a configuration model: the hidden
edges' own degree sequence, rewired. A model that has learned fame scores the rewired pairs as highly
as the real ones and fails. Both easy numbers are still computed and reported beside the hard one, in
`numbers`, because the gap between them is a quantity a reader deserves rather than a detail.

Edges are hidden by orthogroup, never at random, so a hidden edge cannot be recovered through a
paralog that stayed visible; and the layer being predicted is removed from its own features, which is
the rule `graphspace` enforces in one place and these two strategies inherit.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import graphspace as G
from . import strategies as S
from .strategies import Param, Strategy, StrategyResult, register

#: The eighth family, reused rather than invented: 31 and 32 combine strategies, 33 and 34 combine
#: the layers. A new family would have been a ninth heading for the same idea.
COMBINE = "Combine strategies"

EXCLUDE = Param("exclude", "category", "Label being held out",
                "A label this space must stay blind to, if the neighbourhoods are going to be used "
                "while scoring that label. Its whole closure is applied: the column, anything "
                "restating it, the experiment that produced it, and any edge layer built from it.",
                None, optional=True)
KNN = Param("knn", "int", "Nearest genes per gene",
            "How many nearest genes in measurement space become candidate neighbours for each gene, "
            "on top of its partners in every layer. This is what lets the space speak about a pair "
            "no network touches at all; zero restricts it to pairs some layer already links.",
            G.KNN, lo=0, hi=100)
PER_LAYER = Param("per_layer", "int", "Neighbours kept per layer",
                  "How many partners per gene are kept from each layer, strongest weight first. A "
                  "dense layer such as the shared-compartment one would otherwise supply more than "
                  "a hundred thousand pairs that all make the same claim.",
                  G.PER_LAYER, lo=5, hi=200, step=5)


def _layer_default(ctx):
    """The layer the self-test hides: co-expression where it exists, and the reason is measured.

    Crosslinks would be the obvious default -- a missing crosslink is the concrete claim a user acts
    on -- and they are the one layer where this test's null is not a fair contrast. Crosslinked
    complexes are cliques, so a degree-preserving rewiring of the crosslink edges lands mostly inside
    the same complexes, and a pair inside a complex that no layer records looks to every other source
    exactly like one that does. Measured on the shipped table: the model reaches AUROC 0.875 on hidden
    crosslinks against degree-matched non-pairs and the configuration-model null reaches 0.832, so the
    verdict is a FAIL by +0.043 where 0.05 was needed. The other four layers pass, co-expression by
    +0.296 on 3,023 hidden edges. So the default is co-expression -- enough hidden edges for a stable
    verdict against a null that is a real contrast -- and the crosslink layer is one setting away,
    with its failure written down in `docs/graphspace.md` rather than avoided.
    """
    try:
        targets = G.Evidence(ctx).targets()
    except Exception:                       # a table with no layers at all; the runner says so
        return None
    for layer in ("coexpression", "cofitness", "struct", "cotranslation", "xlms"):
        if layer in targets:
            return layer
    return targets[0] if targets else None


LAYER = Param("layer", "layer", "Layer to hide in the self-test",
              "Which layer's edges are hidden and asked back when the strategy tests itself. It is "
              "removed from its own features, so the question is whether the OTHER evidence finds "
              "its edges; a dense layer gives a stricter test than a sparse one.", _layer_default)


def _space(ctx, p, models=G.MODELS) -> G.GraphSpace:
    """The integrated space at this strategy's settings, with the honest evaluation attached."""
    return G.space(ctx, exclude=p.get("exclude") or None, per_layer=int(p["per_layer"]),
                   knn=int(p["knn"]), models=models)


def _nulls_table(report: G.EvaluationReport) -> pd.DataFrame:
    """The three-null comparison, with the fame gap as its own row rather than a footnote."""
    table = report.frame()
    gaps = pd.DataFrame([
        {"model": "logistic", "null": "fame gap (random - degree_matched)",
         "auroc": report.fame_gap},
        {"model": "fame (the usual recipe, diagnostic, never ranked)",
         "null": "fame gap (random - degree_matched)", "auroc": report.fame_recipe_gap},
        {"model": "embedding - logistic", "null": "degree_matched",
         "auroc": report.learned_gain}])
    return pd.concat([table, gaps], ignore_index=True)


# =========================================================================== 33 · the space
def _space_run(ctx, p):
    sp = _space(ctx, p)
    k = int(p["k"])
    wanted, missing = ctx.resolve_genes(p.get("genes") or "")
    if len(wanted):
        focus = wanted
        chosen = f"the {len(focus):,} genes given"
    else:
        # No list: the genes at the top of the space, so the strategy always shows a neighbourhood
        # rather than an empty table and a shrug.
        top = sp.to_frame(top=int(p["top"]))
        focus = np.unique([sp.position(g) for g in list(top["gene_a"]) + list(top["gene_b"])])[:25]
        chosen = "the genes holding the strongest edges in the space, since no gene list was given"
    rows = []
    for i in focus[:200]:
        table = sp.neighbours(int(i), k).rename(
            columns={"gene_id": "neighbour", "product": "neighbour_product"})
        table.insert(0, "gene_id", ctx.gene_ids[int(i)])
        table.insert(1, "product", ctx.product([int(i)])[0])
        rows.append(table)
    neighbourhoods = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    report = sp.report
    summary = (
        f"{sp.summary()} Neighbourhoods are listed for {chosen}"
        + (f"; {len(missing)} accessions were not found ({', '.join(missing[:5])})."
           if missing else ".")
        + f" Held out by orthogroup and scored against degree-matched non-pairs, the space reaches "
          f"AUROC {report.auroc('logistic', 'degree_matched'):.3f}; against random non-pairs "
          f"{report.auroc('logistic', 'random'):.3f}. A neighbour marked measured is an edge some "
          f"experiment recorded; one marked inferred is this space's claim and nothing more.")
    return StrategyResult("neighbour_space", summary,
                          {"neighbourhoods": neighbourhoods,
                           "the space": sp.to_frame(top=int(p["top"])),
                           "how each source is weighted": report.weights(),
                           "the three nulls": _nulls_table(report),
                           "per held-out layer": report.per_layer()},
                          genes=[str(g) for g in ctx.gene_ids[focus][:200]],
                          numbers={"candidate_pairs": len(sp.pairs),
                                   "inferred_pairs": int((~sp.measured).sum()),
                                   "auroc_degree_matched":
                                       report.auroc("logistic", "degree_matched"),
                                   "auroc_random": report.auroc("logistic", "random"),
                                   "fame_gap": report.fame_gap})


def _space_test(ctx, p):
    return G.hidden_edge_test(ctx, "neighbour_space", layer=p["layer"],
                              exclude=p.get("exclude") or None, model="logistic")


register(Strategy(
    key="neighbour_space", number=33, family=COMBINE,
    title="Put every layer into one space and read a gene's neighbourhood",
    method="logistic edge model",
    question="Which genes are the nearest neighbours of this one when every permitted network and "
             "the whole measurement table are combined into a single graph, and what evidence puts "
             "each of them there?",
    tooltip="Builds one graph over the genes from every permitted edge layer plus nearest "
            "neighbours in measurement space, learns what each source is worth from held-out edges, "
            "and lists a gene's neighbours with the probability and the per-source evidence behind "
            "each; measured edges stay marked as measured.",
    explanation=(
        "The thirteen edge layers are deliberately never merged, because each answers a different "
        "question and they disagree. That is right for reading a layer and wrong for the question a "
        "user actually asks: which genes are near this one? Answering it needs all of them at once. "
        "This strategy builds that graph. Candidate pairs are every gene's partners in every "
        "permitted layer plus its nearest genes in measurement space, so a pair no network touches "
        "can still be a neighbour. For each pair the sources are each layer's weight as a percentile "
        "among that layer's own edges, an indicator that the layer does not record the pair at all, "
        "the cosine similarity of the permitted measurements, and an Adamic-Adar shared-partner "
        "count that discounts partners shared through a hub.\n\n"
        "Each source's weight is learned from held-out edges, one fold per layer, and in a layer's "
        "own fold that layer is not among its features -- so a source that can only predict itself "
        "ends at weight zero. Degree is not a source at all: it is the strongest predictor of an "
        "edge and the least interesting, and a model given it learns which genes are famous. The "
        "result is one probability per pair with the per-source contributions that produced it, and "
        "every table marks whether a pair is measured or inferred. Nothing here is written back into "
        "a layer: an integrated edge is a claim, and a crosslink is an experiment."),
    walkthrough=(
        "Paste the genes whose neighbourhoods you want, or leave the list empty to see the "
        "strongest neighbourhoods in the space.",
        "If you are going to use this while scoring a label, name that label under 'Label being "
        "held out' so its whole closure is removed first.",
        "Press Test: the AUROC of hidden edges against DEGREE-MATCHED non-pairs, with the easy "
        "random-null number reported beside it so you can see the difference.",
        "Press Run and read 'neighbourhoods': 'measured_in' names the layers that record the pair "
        "and is empty for an inferred neighbour; 'top_evidence' says which sources moved it.",
        "Read 'the three nulls' before quoting any number from this strategy, and 'how each source "
        "is weighted' to see which kinds of evidence the table actually rests on."),
    test_description=(
        "A layer's edges are hidden by orthogroup -- whole groups at a time, so no hidden edge "
        "survives through a visible paralog -- and the layer is removed from its own features. "
        "Metric: AUROC of the hidden edges against DEGREE-MATCHED non-pairs, one per hidden edge, "
        "matched on connectivity at both ends. Null: 10 configuration-model rewirings of the hidden "
        "edges, which keep their degree sequence and destroy their topology, so a model reading "
        "fame cannot beat it. Pass: above the null's 95th percentile by 0.05. The AUROC against "
        "random non-pairs and the gap between the two are reported as numbers, not as the verdict."),
    params=(Param("genes", "genes", "Genes to read",
                  "Gene accessions whose neighbourhoods to list, one per line or separated by "
                  "commas or spaces. Leave empty to list the neighbourhoods around the strongest "
                  "edges in the space instead. Anything not found is reported by name.", ""),
            EXCLUDE, LAYER,
            Param("k", "int", "Neighbours per gene",
                  "How many neighbours to list for each gene given. The space ranks every candidate "
                  "pair; this only decides how far down each gene's list is printed.", 10, lo=1,
                  hi=200),
            Param("top", "int", "Pairs to list",
                  "How many of the strongest pairs in the whole space to list in 'the space', and "
                  "how many to draw the fallback genes from when no gene list was given.", 500,
                  lo=10, hi=50000, step=50),
            KNN, PER_LAYER),
    runner=_space_run, tester=_space_test, cost="minutes",
    needs=("at least one edge layer",)))


# =========================================================================== 34 · the training
def _training_run(ctx, p):
    sp = _space(ctx, p)
    report = sp.report
    wanted = str(p["model"])
    earned = report.learned_earns_its_place
    used = "logistic"
    if wanted == "embedding" and earned:
        used = "embedding"
    gaps = sp.gaps(int(p["top"]) * 4)
    gaps = gaps[gaps["probability"] >= float(p["min_probability"])].head(int(p["top"]))
    gaps = gaps.reset_index(drop=True)
    note = ""
    if wanted == "embedding" and not earned:
        note = (f"The learned embedding was asked for and is NOT used to rank: on this table it adds "
                f"{report.learned_gain:+.3f} AUROC on the degree-matched null, which does not clear "
                f"the bar (0.01, and twice its own spread across the held-out layers). It was still "
                f"fitted and its numbers are in 'the three nulls'. The ranking is the logistic "
                f"baseline, which can say why it ranked each pair where it did. ")
    elif wanted == "embedding":
        note = (f"The learned embedding earns its place on this table ({report.learned_gain:+.3f} "
                f"AUROC on the degree-matched null) and was used for the ranking; its per-pair "
                f"evidence is the baseline's, since an embedding has no per-source attribution. ")
    summary = (
        f"Trained on {len(report.scored)} held-out layers "
        f"({', '.join(f.layer for f in report.scored)}), each blind to itself, over "
        f"{len(sp.layers)} permitted layers and {len(sp.columns)} permitted measurements. "
        f"Held-out edges reach AUROC {report.auroc(used, 'degree_matched'):.3f} against "
        f"degree-matched non-pairs, {report.auroc(used, 'random'):.3f} against random ones "
        f"(fame gap {report.fame_gap:+.3f}; the same recipe trained the usual way -- degree "
        f"as a feature, random non-pairs to train against -- has a gap of "
        f"{report.fame_recipe_gap:+.3f}), Brier "
        f"{report.metric('brier', used, 'degree_matched'):.3f} with a reliability gap of "
        f"{report.metric('reliability_gap', used, 'degree_matched'):.3f}. {note}"
        f"'gaps' lists {len(gaps):,} pairs the combined evidence implies at probability "
        f">= {float(p['min_probability']):.2f} that NO permitted layer records -- each one a "
        f"relationship this data asserts and no experiment here has seen. Every gap's probability "
        f"rests on similarity, shared partners and shared annotation only, since a pair no layer "
        f"records has no layer weight to rest on.")
    return StrategyResult("network_training", summary,
                          {"gaps": gaps,
                           "the three nulls": _nulls_table(report),
                           "per held-out layer": report.per_layer(),
                           "how each source is weighted": report.weights(),
                           "measured edges, ranked": sp.to_frame(top=int(p["top"]),
                                                                 measured=True)},
                          genes=list(dict.fromkeys(list(gaps.get("gene_a", []))
                                                   + list(gaps.get("gene_b", []))))[:200],
                          numbers={"model_used": used, "learned_gain": report.learned_gain,
                                   "learned_earns_its_place": bool(earned),
                                   "auroc_degree_matched": report.auroc(used, "degree_matched"),
                                   "auroc_random": report.auroc(used, "random"),
                                   "fame_gap": report.fame_gap,
                                   "fame_gap_usual_recipe": report.fame_recipe_gap,
                                   "brier": report.metric("brier", used, "degree_matched"),
                                   "gaps": len(gaps)})


def _training_test(ctx, p):
    return G.hidden_edge_test(ctx, "network_training", layer=p["layer"],
                              exclude=p.get("exclude") or None, model=str(p["model"]),
                              fraction=float(p["fraction"]))


register(Strategy(
    key="network_training", number=34, family=COMBINE,
    title="Train on the networks and rank the edges they are missing",
    method="logistic / spectral embedding",
    question="Which pairs of genes does the combined evidence imply although no measured layer "
             "records them, how strong is each claim, and how good is the model that makes it when "
             "it is scored against a degree-matched null rather than a random one?",
    tooltip="Trains one model per layer, each blind to its own layer, on edges held out by "
            "orthogroup; reports every number against three nulls including degree-matched "
            "non-pairs; and ranks the pairs the evidence implies that no layer records, each with "
            "its calibrated probability and the sources that produced it.",
    explanation=(
        "A network is a sample of a biology, and the pairs it never recorded include real ones. "
        "Which ones can be said from everything else: the other layers, the measurements, who is "
        "whose neighbour. This strategy learns that, one fold per layer, and holds two rules "
        "throughout. Edges are held out by ORTHOGROUP, whole groups at a time, because a randomly "
        "held-out edge is recovered through a paralog that stayed visible, which scores the "
        "resemblance of a gene to its own copy. And a layer is never a feature for predicting "
        "itself, which is enforced in one place and reflected in the weights: a source's weight is "
        "the average of what it was worth in every fold, counting as zero the fold where it was the "
        "target.\n\n"
        "Every number is reported against three nulls at once, because the usual one flatters. "
        "Against random non-pairs a link predictor mostly rediscovers that well-connected genes are "
        "well connected -- strategy 16 reports AUROC 0.97 that way. Against DEGREE-MATCHED "
        "non-pairs, each positive paired with a non-pair of similar connectivity at both ends, the "
        "honest number is smaller and it is the one reported here; against a configuration model, "
        "the degree sequence rewired, a model that has only learned fame collapses to chance. The "
        "difference between the first two is printed as the fame gap, and a deliberately "
        "disqualified degree-only model is fitted beside the real one to show what that gap looks "
        "like when fame IS the model. Probabilities are calibrated on held-out pairs and their "
        "reliability is reported, because a gap table invites a reader to act on a number.\n\n"
        "One layer defeats this test on the shipped data and the failure is the useful part. Hidden "
        "crosslinks are recovered at AUROC 0.875 against degree-matched non-pairs -- and a "
        "configuration-model rewiring of those same crosslinks is recovered at 0.832, because a "
        "crosslinked complex is a clique and a rewired pair inside one looks, to every source except "
        "the crosslink layer itself, exactly like a recorded pair. So this space can say which "
        "COMPLEXES the crosslinking missed and cannot say which PAIRS inside them it missed. Read a "
        "predicted crosslink that way.\n\n"
        "A learned spectral embedding is fitted beside the interpretable baseline and offered for "
        "ranking only where it beats it on the degree-matched null by more than 0.01 AUROC and by "
        "more than twice its own spread across the held-out layers. On the shipped Toxoplasma table "
        "it does not, and the strategy says so and ranks with the baseline: a gap table without "
        "per-source attribution is worth less than a slightly worse one that can be read."),
    walkthrough=(
        "Leave the model on the logistic baseline unless you want to see what the learned "
        "embedding does on your table; it is always fitted and always reported either way.",
        "If you are going to use the result while scoring a label, name that label under 'Label "
        "being held out' so its closure is removed before anything is computed.",
        "Press Test: hidden edges against degree-matched non-pairs, with a configuration-model bar. "
        "Read the 'auroc_random' number beside it to see how much the easy null would flatter.",
        "Press Run and read 'the three nulls' FIRST. If the degree-matched AUROC is near 0.5, the "
        "gaps below it are not worth reading whatever their probabilities say.",
        "Then read 'gaps': each row is a pair no layer records, with its probability and the "
        "sources that produced it. Take the top ones to strategy 13 or 16 for a physical reading.",
        "Check 'how each source is weighted' to see which evidence the ranking rests on, and "
        "whether a source you do not trust on this table is carrying it."),
    test_description=(
        "The chosen layer's edges are hidden by orthogroup and the layer is removed from its own "
        "features. Metric: AUROC of the hidden edges against DEGREE-MATCHED non-pairs. Null: 10 "
        "configuration-model rewirings of the hidden edges, which preserve their degree sequence, "
        "so fame alone cannot clear the bar. Pass: above the null's 95th percentile by 0.05. The "
        "random-null AUROC, the fame gap, precision@k, the Brier score and the reliability gap are "
        "all reported as numbers beside the verdict."),
    params=(EXCLUDE, LAYER,
            Param("model", "choice", "Model",
                  "The logistic baseline on the pair features is interpretable and always fitted. "
                  "The embedding adds spectral node vectors and is used for ranking only if it "
                  "beats the baseline on the degree-matched null; otherwise it is reported and the "
                  "baseline ranks.", "logistic", choices=("logistic", "embedding")),
            Param("top", "int", "Gaps to list",
                  "How many of the highest-probability pairs that no layer records to list. Every "
                  "candidate pair is scored; this only limits what is shown and saved.", 300,
                  lo=10, hi=20000, step=50),
            Param("min_probability", "float", "Minimum probability",
                  "Gaps below this calibrated probability are not listed. The probability is "
                  "against a degree-matched non-pair, not an absolute posterior, so treat it as a "
                  "ranking with a scale rather than as a chance of being true.", 0.0, lo=0.0,
                  hi=1.0, step=0.05),
            Param("fraction", "float", "Share of orthogroups held out",
                  "How much of the graph is hidden for scoring, by orthogroup. A larger share is a "
                  "harder test on less training data; a smaller one leaves too few hidden edges in "
                  "a sparse layer to score at all.", G.TEST_FRACTION, lo=0.1, hi=0.5, step=0.05),
            KNN, PER_LAYER),
    runner=_training_run, tester=_training_test, cost="minutes",
    needs=("at least one edge layer",)))
