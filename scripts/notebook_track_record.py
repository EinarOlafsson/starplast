"""Write notebooks/track_record_2026_10_03.ipynb: what the shipped hold-out track record says.

    python scripts/notebook_track_record.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook("The hold-out track record")
nb.md("# The hold-out track record",
      "Every labelled gene was hidden once and each strategy that calls labels was asked what it is. "
      "This notebook reads the record that ships with Starplast and states what it says. It was "
      "built by `scripts/build_track_record.py` (`starplast/track_record.py`), design in "
      "`instructions/open/62_holdout_track_record.md`.",
      "Three rules make it trustworthy rather than merely present: the five folds never split an "
      "orthogroup (a paralogue cannot give away its sibling's label); a strategy that says nothing "
      "is counted apart from one that is wrong; and no percentage is printed from fewer than five "
      "answered genes.")
nb.code("import pandas as pd",
        "from starplast import track_record as T",
        "led = pd.read_parquet(T.shipped_path())",
        "led.groupby(['organism', 'mode'], observed=True).size().rename('rows')")
nb.md("## Every recorded label",
      "Each organism's biological labels (`track_record.labels`): the declared targets, and for "
      "Toxoplasma the cell-cycle phase and four specific screen phenotypes. Provenance columns are "
      "left out on purpose. Per label: the best strategy on mean per-class recall (small classes "
      "count equally), its overall accuracy, and what always naming the commonest class scores.")
nb.code("rows = []",
        "for (org, target), part in led[led['mode'] == 'together'].groupby(['organism', 'target'], observed=True):",
        "    every = part.drop_duplicates('gene')['truth'].astype(str)",
        "    best = max(((T._baselines(p), k) for k, p in part.groupby(part['strategy'].astype(str))",
        "                if T._baselines(p) is not None), key=lambda t: t[0]['balanced'])",
        "    rows.append({'organism': org, 'label': target, 'genes': len(every), 'classes': every.nunique(),",
        "                 'best strategy': best[1], 'per class': best[0]['balanced'],",
        "                 'chance': 1 / every.nunique(), 'overall': best[0]['accuracy'],",
        "                 'commonest': every.value_counts(normalize=True).iloc[0]})",
        "pd.DataFrame(rows).round(2)")
nb.md("## Coverage",
      "`together` rows are the fold hold-outs: one per (strategy, labelled gene). `set` rows hide "
      "several genes at once. Each fold-held gene must appear once per strategy -- no more, no less.")
nb.code("folds = led[led['mode'] == 'together']",
        "per = folds.groupby(['organism', 'target', 'strategy'], observed=True)['gene'].agg(['size', 'nunique'])",
        "assert (per['size'] == per['nunique']).all(), 'a gene was held out twice'",
        "per.rename(columns={'size': 'held out', 'nunique': 'distinct genes'})")
nb.md("The strategy and class levels below are each organism's main label (`T.default_target`); the "
      "table above covers the rest.")
nb.code("main = pd.concat([folds[(folds['organism'] == o) & (folds['target'].astype(str) == T.default_target(o, led))]",
        "                  for o in folds['organism'].astype(str).unique()])",
        "main.groupby('organism', observed=True)['target'].first()")
nb.md("## Strategy level",
      "Right over answered, with the Wilson 95% interval. `answered` below the number held out means "
      "the strategy abstained on the rest -- shown, not hidden in the rate.")
nb.code("s = T.summary(main, 'strategy')",
        "s[['organism', 'strategy', 'answered', 'right', 'wrong', 'abstained', 'rate', 'rate_low', 'rate_high']]"
        ".sort_values(['organism', 'rate'], ascending=[True, False]).round(3)")
nb.md("## Class level: which distinctions the data carries",
      "For each class, the best strategy's rate. A class every strategy misses is a statement about "
      "the measurements, not about any one strategy.")
nb.code("c = T.summary(main[main.organism == 'Tg'], 'class')",
        "c = c[c['enough'].astype(bool)]",
        "best = c.sort_values('rate', ascending=False).drop_duplicates('truth')",
        "best[['truth', 'strategy', 'right', 'answered', 'rate', 'confused_with']].round(3)")
nb.md("## Hidden together: does it still work when many genes are unknown?",
      "Random sets of 1 to 500 genes hidden at once. A strategy that only works when a gene's "
      "neighbours are labelled falls as the set grows; a flat curve means the evidence reaches "
      "the genes themselves.")
nb.code("sets = led[(led['mode'] == 'set') & led['set_name'].astype(str).str.startswith('random')]",
        "curve = sets.assign(answered=~sets['abstained'].astype(bool))",
        "curve = curve[curve['answered']].groupby(['organism', 'target', 'strategy', 'set_size'], observed=True)['correct']"
        ".mean().unstack('set_size')",
        "curve.round(2)")
nb.md("## Whole classes hidden",
      "Every gene of one class hidden at once. `right` is zero by construction: with every example "
      "of a label hidden, no strategy that learns from examples can name it. What is measured is "
      "`together` -- the share given one and the same label -- i.e. whether the strategy still "
      "recognises them as one thing, and `placed_at`, where it put them. Only sets with five or "
      "more answered members are ranked.")
nb.code("cls = led[(led['mode'] == 'set') & ~led['set_name'].astype(str).str.startswith('random')]",
        "w = T.set_summary(cls)",
        "assert (w['right'] == 0).all(), 'a hidden label was predicted: the class was not fully hidden'",
        "w = w[w['answered'] >= 5]",
        "w[['organism', 'strategy', 'set_name', 'set_size', 'answered', 'together', 'placed_at']]"
        ".sort_values('together', ascending=False).head(20).round(2)")
nb.md("## Your own list",
      "`T.my_list(ctx, genes, target)` hides a pasted list together and asks every strategy. "
      "Two real labelled genes, as an example:")
nb.code("from starplast import strategies as S",
        "ctx = S.Context.shipped('Tg')",
        "rows = T.my_list(ctx, ['TGME49_294550', 'TGME49_244470'], 'compartment', "
        "strategies=['feature_knn', 'physical_partners', 'layer_vote'])",
        "rows[['strategy', 'gene_id', 'truth', 'prediction', 'correct', 'abstained']]")
nb.write(os.path.join(ROOT, "notebooks", "track_record_2026_10_03.ipynb"))
