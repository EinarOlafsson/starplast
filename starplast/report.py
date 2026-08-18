#!/usr/bin/env python3
"""The figures a recipe run would go into a paper with.

Instruction 46. The application already lets a result be read on screen; this writes the same run out
as a document somebody could put in a manuscript -- vector output, one page per map, and every number
on it traceable to a rebuildable run.

**One page per UMAP, and the alternatives are included on purpose.** A run tunes several maps and
keeps the best; a report that showed only the winner would hide what the choice was between, which is
exactly the thing a reviewer asks about. Each page reads left to right in the order the question is
actually asked: what the map's settings were, what the clustering's settings were, the map coloured by
cluster, by the label it had to recover, and by the independent control -- then, underneath, how well
each of those labels actually maps onto those clusters.

**The three-way comparison is the point.** Cluster against label says whether the structure found the
biology; cluster against validation label says whether it found the biology it was ASKED about rather
than something else that happens to separate genes. Those two panels sit above their own quality
graphs so the eye can move between the picture and the number.

**A 3D embedding drawn on paper is a projection, and the page says so.** The scatter panels are
components 1 and 2 of a three-component UMAP, labelled as such, because a reader who assumes they are
seeing the whole map will read two clusters as touching when they are separated in the third.
"""
from __future__ import annotations

import os
from datetime import datetime

import numpy as np
import pandas as pd

#: How many clusters the label-by-cluster heatmap shows. The full matrix is labels x clusters and a
#: tuned map has over a hundred of the latter, which prints as a grey smear. The largest clusters are
#: the ones a claim can be made about, so those are the ones drawn -- and the caption says how many
#: were left out rather than letting the picture imply there were none.
HEATMAP_CLUSTERS = 24

#: Distinct labels drawn in colour on a scatter before the rest become "other". Beyond this the eye
#: cannot tell the colours apart and the legend is taller than the figure.
SCATTER_LABELS = 8


def _palette(n: int):
    from matplotlib import colormaps
    base = list(colormaps["tab20"].colors) + list(colormaps["tab20b"].colors)
    return [base[i % len(base)] for i in range(max(n, 1))]


def _limits(coords, span: float = 99.0, pad: float = 0.04):
    """Axis limits from the bulk of the cloud rather than from its outliers.

    Measured on a real run: a UMAP of this table put a handful of genes at x = 10 while the other
    8,100 sat between -52 and -46, and matplotlib's autoscale duly drew all the structure as a
    four-pixel smear against an empty page. The percentile range is what a reader needs to see; the
    outliers are still plotted, just outside the frame.
    """
    coords = np.asarray(coords)
    if not len(coords):
        return (-1.0, 1.0), (-1.0, 1.0)
    lo = np.percentile(coords[:, :2], (100 - span) / 2, axis=0)
    hi = np.percentile(coords[:, :2], 100 - (100 - span) / 2, axis=0)
    width = np.maximum(hi - lo, 1e-6)
    return ((lo[0] - pad * width[0], hi[0] + pad * width[0]),
            (lo[1] - pad * width[1], hi[1] + pad * width[1]))


def _scatter(ax, coords, values, title: str, noise_label=None):
    """One embedding, coloured by a categorical column, with absence drawn as absence.

    Unlabelled genes are drawn first, in pale grey, rather than omitted: a cluster that is mostly
    unmeasured looks identical to a small cluster if the unlabelled genes are simply not there, and
    which is which decides whether the cluster supports a claim.
    """
    coords = np.asarray(coords)
    values = pd.Series(values).astype("object").where(pd.notna(pd.Series(values)), None).to_numpy()
    # dtype=bool explicitly: an empty list becomes a float64 array, and `~` on that raises rather
    # than returning an empty mask. A map with no genes is a real case -- every gene can fail the
    # missing-value policy -- and it must draw an empty panel, not crash the document.
    known = np.array([v is not None and str(v) != str(noise_label) for v in values], dtype=bool)
    ax.scatter(coords[~known, 0], coords[~known, 1], s=2.5, c="#dedede", linewidths=0,
               label=f"unlabelled ({int((~known).sum())})", rasterized=True)
    counts = pd.Series([str(v) for v in values[known]]).value_counts() if known.any() \
        else pd.Series(dtype=int)
    top = list(counts.index[:SCATTER_LABELS])
    colors = _palette(len(top))
    for name, color in zip(top, colors):
        m = known & np.array([str(v) == name for v in values])
        ax.scatter(coords[m, 0], coords[m, 1], s=5.0, color=color, linewidths=0,
                   label=f"{name} ({int(m.sum())})", rasterized=True)
    rest = known & np.array([str(v) not in top for v in values], dtype=bool)
    if rest.any():
        ax.scatter(coords[rest, 0], coords[rest, 1], s=2.0, c="#7a7a7a", linewidths=0,
                   label=f"other ({int(rest.sum())})", rasterized=True)
    ax.set_title(title, fontsize=8)
    ax.set_xlabel("UMAP 1 of 3", fontsize=6)
    ax.set_ylabel("UMAP 2 of 3", fontsize=6)
    ax.tick_params(labelsize=5)
    xlim, ylim = _limits(coords)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    # Inside the axes, and this is the second attempt. Placed BELOW the panel the legends ran into
    # the row underneath and sat on top of its titles and its summary text -- a figure that hides
    # another figure is worse than one that hides some of its own points. `best` lets matplotlib put
    # it wherever the cloud is thinnest, which on a UMAP is usually a corner.
    ax.legend(fontsize=4, markerscale=3, ncol=1, framealpha=0.88, loc="best", borderpad=0.4)


def _table(ax, rows, title: str):
    """A settings table as a figure panel, because a figure without its settings is not evidence."""
    ax.axis("off")
    ax.set_title(title, fontsize=8, loc="left")
    # No empty-rows guard: every caller passes a fixed list of settings, so an empty table cannot
    # occur. A branch that cannot run is deleted rather than tested -- the alternative is declaring
    # it covered.
    body = [[str(k), str(v)] for k, v in rows]
    t = ax.table(cellText=body, colWidths=[0.52, 0.48], loc="upper left", cellLoc="left")
    t.auto_set_font_size(False)
    t.set_fontsize(6)
    t.scale(1, 1.25)
    for cell in t.get_celld().values():
        cell.set_linewidth(0.3)


def _per_label(ax, recovery: pd.DataFrame, title: str):
    """Precision and recall per label, against the single cluster that recovered it best."""
    ax.set_title(title, fontsize=8)
    if not len(recovery):
        ax.axis("off")
        ax.text(0.5, 0.5, "nothing scoreable", ha="center", va="center", fontsize=7)
        return
    r = recovery.sort_values("f1", ascending=True).tail(14)
    y = np.arange(len(r))
    ax.barh(y - 0.2, r.precision, height=0.38, color="#4c72b0", label="precision")
    ax.barh(y + 0.2, r.recall, height=0.38, color="#dd8452", label="recall")
    ax.set_yticks(y)
    ax.set_yticklabels([str(v)[:26] for v in r.label], fontsize=5)
    ax.set_xlim(0, 1)
    ax.tick_params(labelsize=5)
    ax.legend(fontsize=5, loc="lower right")


def _heatmap(ax, per_cluster: pd.DataFrame, title: str):
    """The full label x cluster matrix, which the best-cluster summary cannot show.

    A label split cleanly across three clusters is a finding about sub-structure; read only through
    its best cluster it looks like a label that half-recovered.
    """
    ax.set_title(title, fontsize=8)
    if not len(per_cluster):
        ax.axis("off")
        ax.text(0.5, 0.5, "nothing scoreable", ha="center", va="center", fontsize=7)
        return
    keep = (per_cluster.groupby("cluster").n_labelled_in_cluster.max()
            .sort_values(ascending=False).head(HEATMAP_CLUSTERS).index)
    m = per_cluster[per_cluster.cluster.isin(keep)].pivot_table(
        index="label", columns="cluster", values="f1", fill_value=0.0)
    im = ax.imshow(m.values, aspect="auto", cmap="magma", vmin=0, vmax=max(m.values.max(), 0.01))
    ax.set_xticks(range(len(m.columns)))
    ax.set_xticklabels(m.columns, fontsize=4, rotation=90)
    ax.set_yticks(range(len(m.index)))
    # Truncated from the LEFT would lose the distinguishing end; these names differ at the start, so
    # they are cut at 16 and the panel keeps its margin instead of running off the page.
    ax.set_yticklabels([str(v)[:16] for v in m.index], fontsize=4)
    ax.set_xlabel(f"cluster (largest {len(m.columns)} of "
                  f"{per_cluster.cluster.nunique()})", fontsize=6)
    ax.figure.colorbar(im, ax=ax, fraction=0.035).ax.tick_params(labelsize=4)


def _agreement(ax, control: pd.DataFrame, title: str):
    """Where the control landed on the primary's clusters, and whether it corroborates them."""
    ax.set_title(title, fontsize=8)
    # The column is `n_control_labelled` here: `run` renames it when it stores the control table, so
    # that a merge onto the inference cannot collide with the primary's own counts. Reading
    # `n_labelled` instead crashed on every real run that had a control, and passed every test built
    # from `dominant_by_cluster` directly -- which is why the fixture now builds it the way `run`
    # does, renames and all.
    count = "n_control_labelled" if "n_control_labelled" in control.columns else "n_labelled"
    scoreable = control[control[count] > 0] if len(control) else control
    if not len(scoreable):
        ax.axis("off")
        ax.text(0.5, 0.5, "the control has no gene in any cluster:\nit cannot corroborate or refute",
                ha="center", va="center", fontsize=7)
        return
    r = scoreable.sort_values("enrichment", ascending=False).head(24)
    colors = ["#55a868" if a else "#c44e52" for a in r.agrees]
    ax.bar(range(len(r)), r.enrichment, color=colors)
    ax.axhline(1.0, color="#333333", lw=0.6, ls="--")
    # Log scale, because enrichment is a ratio and a rare control label reaches several hundred while
    # a common one sits near 2. Linear, the informative end of the axis is one pixel tall.
    if float(r.enrichment.max()) / max(float(r.enrichment.min()), 1e-9) > 20:
        ax.set_yscale("log")
    ax.set_xticks(range(len(r)))
    ax.set_xticklabels(r.cluster, fontsize=4, rotation=90)
    ax.set_xlabel("cluster", fontsize=6)
    ax.set_ylabel("control enrichment", fontsize=6)
    ax.tick_params(labelsize=5)
    ax.set_title(f"{title} — {int(r.agrees.sum())} of {len(r)} shown corroborate", fontsize=8)


def _disagreement(ax, result):
    """The clusters an answer actually rests on, and whether the control backs each one.

    The interesting labels are not the best-recovered ones -- those are already in the panel to the
    left -- but the clusters the inference drew genes from, split by whether the control agrees. A
    run whose named genes come entirely from red bars is a run standing on its primary alone.
    """
    ax.set_title("clusters the answer rests on", fontsize=8)
    inference = result.inference
    if not len(inference):
        ax.axis("off")
        ax.text(0.5, 0.5, "no genes were named,\nso no cluster carries the answer",
                ha="center", va="center", fontsize=7)
        return
    per = (inference.groupby(["cluster", "predicted"])
           .agg(genes=("gene_id", "size"), purity=("cluster_precision", "max"),
                enrichment=("enrichment", "max"),
                agrees=("control_agrees", "max") if "control_agrees" in inference else
                       ("cluster_precision", lambda s: False))
           .reset_index().sort_values("genes", ascending=False).head(14))
    y = np.arange(len(per))
    colors = ["#55a868" if bool(a) else "#c44e52" for a in per.agrees]
    ax.barh(y, per.genes, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{str(p)[:20]} (c{int(c)}, {q:.0%})"
                        for p, c, q in zip(per.predicted, per.cluster, per.purity)], fontsize=5)
    ax.invert_yaxis()
    ax.set_xlabel("genes named  (green = control corroborates)", fontsize=6)
    ax.tick_params(labelsize=5)


def _cover(pdf, result, version: str):
    """The page that makes the rest of the document evidence rather than pictures."""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(11.7, 8.3))
    fig.text(0.06, 0.94, "starplast — recipe report", fontsize=15, weight="bold")
    fig.text(0.06, 0.90, result.recipe.question, fontsize=11, wrap=True)
    closure = result.closure
    removed = "\n".join(f"    {c} — {why}" for c, why in sorted(closure.removed.items())[:12]) \
        or "    nothing the question named had to be removed"
    body = (
        f"axis                {result.recipe.axis or '—'}\n"
        f"holdout             {closure.holdout_column}"
        f"{f' (binned into {result.recipe.holdout_bins})' if result.recipe.holdout_bins else ''}\n"
        f"validation holdout  {closure.control_column or 'none available'}\n"
        f"gate                enrichment >= {result.recipe.min_enrichment:.2f}, "
        f"purity >= {result.recipe.min_precision:.2f}\n"
        f"seed                {result.recipe.seed}\n"
        f"leakage scope       {result.recipe.scope}\n\n"
        f"inputs after closure  {len(closure.blocks)} blocks, {len(closure.columns)} columns\n"
        f"excluded overall      {len(closure.excluded)} columns\n"
        f"removed from the inputs this question named ({len(closure.removed)}):\n{removed}\n\n"
        f"maps in this document {1 + len(result.alternatives)}\n"
        f"genes named           {len(result.inference)}\n\n"
    )
    import textwrap
    body += "expectation\n" + "\n".join(
        "    " + line for line in textwrap.wrap(result.recipe.expectation or "—", 96))
    fig.text(0.06, 0.86, body, fontsize=8, family="monospace", va="top")
    fig.text(0.06, 0.04, f"starplast {version} · generated {datetime.now():%Y-%m-%d %H:%M} · "
                         f"every figure below is from the run described here", fontsize=6)
    pdf.savefig(fig)
    plt.close(fig)


def _map_page(pdf, result, index: int, total: int):
    """One map: its settings, the three colourings, and the quality of each."""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(20, 8.5))
    grid = fig.add_gridspec(2, 5, height_ratios=[1.10, 1.0], hspace=0.42, wspace=0.30,
                            left=0.055, right=0.985, top=0.90, bottom=0.08)
    q = result.quality
    fig.suptitle(f"map {index} of {total} — {result.recipe.question}", fontsize=11, y=0.97)

    _table(fig.add_subplot(grid[0, 0]), [
        ("n_neighbors", result.settings.get("n_neighbors", "—")),
        ("min_dist", result.settings.get("min_dist", "—")),
        ("components", np.asarray(result.coords).shape[1] if len(result.coords) else "—"),
        ("metric", "euclidean"),
        ("seed", result.recipe.seed),
        ("backend", result.summary.get("backend", "—")),
        ("genes mapped", int(np.sum(result.genes))),
    ], "UMAP")
    _table(fig.add_subplot(grid[0, 1]), [
        ("algorithm", "hdbscan"),
        ("selection", "leaf"),
        ("min_cluster_size", result.settings.get("min_cluster_size", "—")),
        ("min_samples", result.settings.get("min_samples", "auto")),
        ("clusters", q.get("clusters", 0)),
        ("clustered", f"{100 * q.get('clustered', 0):.0f}%"),
        ("evenness", f"{q.get('evenness', 0):.2f}"),
        ("degenerate", q.get("why_not") or "no"),
    ], "clustering")

    coords = np.asarray(result.coords)
    cluster_names = pd.Series([f"cluster {v}" if v >= 0 else None for v in result.labels])
    _scatter(fig.add_subplot(grid[0, 2]), coords, cluster_names, "coloured by cluster")
    _scatter(fig.add_subplot(grid[0, 3]), coords, result.truth.to_numpy(),
             f"coloured by holdout: {result.closure.holdout_column}")
    _scatter(fig.add_subplot(grid[0, 4]), coords, result.control_truth.to_numpy()
             if len(result.control_truth) else [None] * len(coords),
             f"coloured by control: {result.closure.control_column or 'none'}")

    _per_label(fig.add_subplot(grid[1, 0]), result.recovery, "holdout: recovery per label")
    _heatmap(fig.add_subplot(grid[1, 1]), result.per_cluster, "holdout: F1 per label x cluster")
    _agreement(fig.add_subplot(grid[1, 2]), result.control, "control on the same clusters")

    _disagreement(fig.add_subplot(grid[1, 3]), result)
    ax = fig.add_subplot(grid[1, 4])
    ax.axis("off")
    ax.set_title("what this map answered", fontsize=8, loc="left")
    named = len(result.inference)
    corr = int(result.inference["control_agrees"].fillna(False).sum()) \
        if "control_agrees" in result.inference else 0
    ax.text(0, 0.95, "\n".join([
        f"labels scored     {result.summary.get('n_labels_scored', 0)}",
        f"mean F1           {result.summary.get('mean_f1', float('nan')):.3f}",
        f"best F1           {result.summary.get('best_f1', float('nan')):.3f}"
        f"  ({result.summary.get('best_label', '—')})",
        f"genes named       {named}",
        f"corroborated      {corr}",
        "",
        (result.inference_note or ("the control corroborates none of these clusters"
                                   if named and not corr else "")),
    ]), fontsize=7, family="monospace", va="top", wrap=True)
    pdf.savefig(fig)
    plt.close(fig)


def recipe_pdf(result, path: str) -> str:
    """Write the run to `path` as a PDF and return the path.

    A refused recipe still gets a document: the refusal and what caused it is the result, and a
    reader looking for the run needs to find the reason rather than an absent file.
    """
    import matplotlib
    matplotlib.use("Agg")                      # the suite has no display, and neither does a server
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from . import __version__
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    maps = [result] + list(result.alternatives)
    with PdfPages(path) as pdf:
        _cover(pdf, result, __version__)
        if not result.ok:
            fig = plt.figure(figsize=(11.7, 8.3))
            fig.text(0.06, 0.8, "no answer", fontsize=15, weight="bold")
            fig.text(0.06, 0.72, result.stopped_because, fontsize=10, wrap=True)
            pdf.savefig(fig)
            plt.close(fig)
            return path
        for i, one in enumerate(maps, 1):
            _map_page(pdf, one, i, len(maps))
    return path
