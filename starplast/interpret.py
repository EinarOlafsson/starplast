#!/usr/bin/env python3
"""Read a run's findings back to the reader, in the order they are worth reading.

A search returns a table. A table of a hundred enriched clusters is not an answer -- it is the same
problem the search was meant to solve, moved one step along, and the step that actually costs a
scientist their afternoon is deciding which four rows matter and what they say. This module does
that step: it ranks findings by how much they are worth acting on and writes each one as the claim
it is, with the numbers that support it and the reason to doubt it.

## What "worth acting on" means here

Four things, multiplied rather than added, because a finding that scores zero on any one of them is
worth nothing whatever the others say.

- **Strength** -- how far the evidence is past correction, and SATURATING. A claim at q = 0.04 and
  one at q = 1e-8 are not the same claim; a claim at q = 1e-8 and one at q = 1e-40 are, because both
  are certain and certainty does not come in degrees a biologist can act on. Without the saturation
  the term runs away with the ranking: measured, a q of 1e-40 about two genes outranked a q of 1e-8
  about twenty-five, which is the opposite of useful. Past `CERTAIN`, only reach and novelty
  separate findings.
- **Reach** -- how many genes it is ABOUT. A perfect enrichment that predicts two genes is a
  confirmation; one that predicts forty is a result.
- **Novelty** -- how unstudied the implicated genes are. Predicting the localisation of a gene with
  two hundred papers is not a prediction. This is the term that turns a statistically strong,
  biologically boring cluster into the low rank it deserves.
- **Independence** -- whether the layer was an input to the map. A circular finding is scored at
  zero and still shown, labelled, because a reader who sees only the clean ones cannot tell whether
  the run was clean.

## What it does not do

It does not decide whether a claim is true, and it says so in every report it writes. Enrichment is
evidence that genes group; whether they group for the reason the label suggests is a question for an
experiment. The wording throughout is "predicted" and "consistent with", never "shows" -- a report
that overstates once is a report nobody can use unread.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

#: Where a claim stops being about "a few genes" and starts being about a set worth a figure.
GOOD_REACH = 12.0

#: Where a q-value stops being evidence and starts being arithmetic. Ten to the minus eight is far
#: past any threshold anyone would argue about, and everything beyond it is treated as equally
#: certain.
CERTAIN = 8.0

#: Publications per gene above which a prediction is not news. Not a cliff: the term is a smooth
#: decay, and a well-studied gene in an otherwise novel cluster does not sink the whole finding.
STUDIED = 8.0


def novelty(nodes: pd.DataFrame, gene_ids) -> float:
    """How unstudied the implicated genes are, in [0, 1]. 1 is a set nobody has published on.

    Measured on publications per gene rather than on whether the label is missing. A gene can be
    unlocalised and still be one of the twenty everybody works on, and predicting its compartment is
    a smaller contribution than predicting one for a gene with no literature at all.
    """
    if not len(gene_ids) or "n_publications" not in nodes.columns:
        return 0.5                            # no literature layer loaded: neither reward nor punish
    key = nodes["gene_id"].astype(str) if "gene_id" in nodes.columns else nodes.index.astype(str)
    hits = nodes.loc[key.isin({str(g) for g in gene_ids}), "n_publications"]
    if hits.empty:
        return 0.5
    papers = pd.to_numeric(hits, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    return float(np.mean(np.exp(-papers / STUDIED)))


def interest(findings: pd.DataFrame, nodes: pd.DataFrame) -> pd.DataFrame:
    """Score and sort findings by how much they are worth reading. Adds the four terms as columns.

    Kept as columns rather than folded into one number and thrown away: a reader who disagrees with
    the weighting -- and a reader who cares about one gene family will disagree -- can re-sort on the
    term they care about instead of arguing with a ranking they cannot see inside.
    """
    if findings is None or findings.empty:
        return pd.DataFrame()
    out = findings.copy()
    q = np.clip(pd.to_numeric(out.get("q", 1.0), errors="coerce").fillna(1.0).to_numpy(), 1e-300, 1)
    out["strength"] = np.clip(-np.log10(q) / CERTAIN, 0.0, 1.0)
    reach = pd.to_numeric(out.get("n_predicted", 0), errors="coerce").fillna(0).to_numpy(float)
    out["reach"] = np.log1p(reach) / np.log1p(GOOD_REACH)
    out["novelty"] = [novelty(nodes, g if isinstance(g, (list, tuple, np.ndarray)) else [])
                      for g in out.get("genes", [[]] * len(out))]
    independent = ~out.get("circular", pd.Series(False, index=out.index)).astype(bool)
    out["interest"] = out.strength * out.reach * out.novelty * independent.astype(float)
    return out.sort_values(["interest", "strength"], ascending=False).reset_index(drop=True)


def _names(nodes: pd.DataFrame, gene_ids, limit: int = 6) -> str:
    """Gene ids with their product names where there is one, as a phrase."""
    if "gene_id" not in nodes.columns:
        return ", ".join(str(g) for g in list(gene_ids)[:limit])
    look = nodes.set_index(nodes["gene_id"].astype(str))
    out = []
    for g in list(gene_ids)[:limit]:
        product = str(look["product"].get(str(g), "")) if "product" in look.columns else ""
        product = "" if product.lower() in ("", "nan", "hypothetical protein") else product
        out.append(f"{g} ({product})" if product else str(g))
    more = len(gene_ids) - len(out)
    return ", ".join(out) + (f", and {more} more" if more > 0 else "")


def sentence(row, nodes: pd.DataFrame) -> str:
    """One finding, written as the claim it is."""
    q = float(row.get("q", 1.0))
    stat = f"q = {q:.1e}" if q < 0.01 else f"q = {q:.3f}"
    cluster = int(row.get("cluster", -1))
    genes = row.get("genes") or []

    if row.get("kind") == "conjunction":
        return (
            f"**Cluster {cluster} is {row['purity']:.0%} {row['category']} by {row['layer']} AND "
            f"{row['other_category']} by {row['other']}** ({row['joint_lift']:.1f}x joint "
            f"enrichment, {row['interaction_ratio']:.1f}x beyond either margin alone, {stat}). "
            f"The crossed structure is a hypothesis generated from two measurements, not proof "
            f"that either state causes the other. Unmeasured members: {_names(nodes, genes)}.")

    if row.get("kind") == "guilt" and row.get("layer_kind") == "discrete":
        return (
            f"**Cluster {cluster}: {row['n_hits']} of {row['n_known']} labelled genes are "
            f"{row['category']}** ({row['purity']:.0%} of the cluster against "
            f"{row['background']:.1%} across the labelled proteome, {row['lift']:.1f}x enrichment, "
            f"{stat}). {int(row['n_predicted'])} genes in the same cluster have no "
            f"{row['layer']} label: {_names(nodes, genes)}. Guilt by association predicts they are "
            f"{row['category']} as well.")

    if row.get("kind") == "guilt":
        direction = "lower" if row.get("effect", 0) < 0 else "higher"
        return (
            f"**Cluster {cluster}: {row['layer']} averages {row['mean_in']:.2f} here against "
            f"{row['mean_out']:.2f} elsewhere** ({direction}, d = {row['effect']:.1f}, {stat}, "
            f"{row['n_known']} measured). {int(row['n_predicted'])} genes in this cluster were "
            f"never measured: {_names(nodes, genes)}. They are predicted to sit at the same "
            f"{direction} end.")

    if row.get("other_kind") == "continuous":
        return (
            f"**Cluster {cluster} is {row['purity']:.0%} {row['category']} by {row['layer']}, and "
            f"its {row['other']} splits in two** ({row['mean_low']:.2f} in {int(row['n_low'])} "
            f"genes against {row['mean_high']:.2f} in {int(row['n_high'])}, a gap of "
            f"{row['gap']:.1f} pooled SD, {stat}). One category, two behaviours -- a subdivision "
            f"inside {row['category']} that {row['layer']} alone does not show. The lower group: "
            f"{_names(nodes, genes)}.")

    groups = row.get("groups") or []
    counts = row.get("counts") or []
    split = ", ".join(f"{g} ({c})" for g, c in zip(groups, counts))
    return (
        f"**Cluster {cluster} is {row['purity']:.0%} {row['category']} by {row['layer']}, but its "
        f"{row['other']} disagrees**: {split} ({stat}). Genes that share one label and not the "
        f"other are where a category is hiding a distinction: {_names(nodes, genes)}.")


def caveats(row, nodes: pd.DataFrame, table: pd.DataFrame = None) -> list:
    """Everything a reader should know before believing a finding. Never empty by design."""
    out = []
    if table is not None and not table.empty and row.get("kind") == "guilt":
        # One cluster, several categories. Seen immediately on the real map: cluster 29 came back
        # 30% ER and 26% golgi, and the report offered the same 113 unlabelled genes to both. The
        # predictions are ALTERNATIVES -- a gene is in one compartment -- and reading them as a list
        # of separate findings quietly doubles the apparent yield of a cluster that is really one
        # observation about the secretory pathway.
        siblings = table[(table.get("cluster") == row.get("cluster"))
                         & (table.get("layer") == row.get("layer"))]
        if len(siblings) > 1:
            others = [str(c) for c in siblings.category if str(c) != str(row.get("category"))]
            out.append(f"this cluster also comes back enriched for {', '.join(others[:4])} -- "
                       f"the predictions are alternatives for the same genes, not additional ones")
    if row.get("circular"):
        out.append(f"CIRCULAR: {row['layer']} was among the measurements this map was built from, "
                   f"so this enrichment is arithmetic rather than evidence.")
    if float(row.get("novelty", 1.0)) < 0.35:
        out.append("the implicated genes are already well published, so the prediction adds little")
    if int(row.get("n_known", 0)) < 20:
        out.append(f"only {int(row.get('n_known', 0))} genes carry the measurement here, so the "
                   f"effect is estimated from few points")
    if float(row.get("q", 1.0)) > 0.01:
        out.append("close to the correction threshold: expect some findings at this level to be "
                   "chance")
    out.append("clustering is evidence that these genes behave alike, not that they do the same "
               "job -- the claim is a hypothesis to test, not a result")
    return out


def report(findings: pd.DataFrame, nodes: pd.DataFrame, top: int = 8, alpha: float = 0.05) -> str:
    """The whole reading, as markdown: what the run found, in order, with the reasons to doubt it.

    Written to be pasted into a lab notebook and understood a month later, which is why every claim
    carries its numbers inline rather than referring to a row of a table that will not be there.
    """
    if findings is None or findings.empty:
        return ("## Nothing to report\n\nNo cluster in this configuration carries a claim that "
                "survives correction. That is a real answer about this map rather than a failure "
                "of the search: try a different set of measurements, or a clustering with smaller "
                "minimum cluster size, before concluding the structure is not there.")
    ranked = interest(findings, nodes)
    live = ranked[(ranked.q <= alpha) & (~ranked.get("circular", False).astype(bool))]
    dead = len(ranked) - len(live)
    lines = ["## What this run found", ""]
    lines.append(f"{len(live)} claims survive correction at q < {alpha} across "
                 f"{len(ranked)} tested"
                 + (f"; {dead} more are either uncorrected noise or circular." if dead else "."))
    kinds = live["kind"].value_counts().to_dict() if len(live) else {}
    if kinds:
        lines.append("")
        lines.append("  ".join(
            f"**{ {'guilt': 'guilt by association', 'disagreement': 'layer disagreement', 'conjunction': 'crossed factors'}.get(k, k)}**: {v}"
            for k, v in kinds.items()))
    for i, (_, row) in enumerate(live.head(top).iterrows(), 1):
        lines += ["", f"### {i}. {sentence(row, nodes)}", ""]
        lines.append(f"*interest {row.interest:.2f} = strength {row.strength:.2f} x reach "
                     f"{row.reach:.2f} x novelty {row.novelty:.2f}*")
        lines += [f"- {c}" for c in caveats(row, nodes, ranked)]
    circular = ranked[ranked.get("circular", pd.Series(False, index=ranked.index)).astype(bool)]
    if len(circular):
        lines += ["", "### Excluded as circular", "",
                  f"{len(circular)} findings concern a layer that helped build this map "
                  f"({', '.join(sorted(set(circular.layer.astype(str))))}). They are shown here so "
                  f"the run can be judged, and they are not evidence about biology."]
    lines += ["", "---", "", "Every claim above is a hypothesis about genes that behave alike. "
              "Clustering cannot distinguish a shared function from a shared route, a shared "
              "assay, or a shared reason for being measured at all."]
    return "\n".join(lines)
