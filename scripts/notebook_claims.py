"""Write notebooks/claims_2026_10_04.ipynb: the measurements behind every shipped claim.

    python scripts/notebook_claims.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook("Claims: generated, then tested independently")
nb.md("# Claims: generated with a measured certainty, then tested independently",
      "The goal (2026-10-04): generate new knowledge from existing data, and test it by means separate "
      "from the inference that generated it -- with every number measured on genes whose answer was "
      "hidden. Design: `instructions/open/63_claims_generate_and_verify.md`; code: "
      "`starplast/claims.py`; build: `scripts/build_claims.py`.",
      "Everything below reads the shipped track record (`starplast/data/track_record.parquet`), where "
      "every labelled gene was hidden once in folds that never split an orthogroup.")
nb.code("import numpy as np, pandas as pd",
        "from starplast import claims as C, strategies as S, track_record as T",
        "tg = T.shipped('Tg'); tg = tg[tg['target'].astype(str) == 'compartment']",
        "len(tg)")
nb.md("## 1. Leakage between generator and verifier, measured",
      "\"Information leakage needs to be quantified not guessed\": a signal peptide is location "
      "information though it shares no dataset with LOPIT. The shared-mistake ratio: when the generator "
      "is wrong, how often the verifier makes the SAME wrong call, over how often it would if the two "
      "were independent given the truth (from the verifier's own confusion matrix). 1.0 = independent. "
      f"A verifier counts as an independent test when the upper end of the 95% bootstrap interval is "
      f"below `C.INDEPENDENT_BELOW`.")
nb.code("pairs = [(g, v) for g in ('feature_knn', 'supervised_classifier', 'random_forest', 'stacking')",
        "         for v in ('physical_partners', 'structural_homology', 'layer_propagation', 'random_forest', 'feature_knn') if g != v]",
        "m = pd.DataFrame([C.shared_mistakes(tg, g, v) for g, v in pairs])",
        "m[['generator', 'verifier', 'generator_wrong', 'observed', 'expected', 'ratio', 'ratio_low', 'ratio_high']].round(2)")
nb.md("Stacking absorbs every kind of evidence, so nothing can test it independently -- which is why "
      "recipes are chosen as the generator that independent evidence CAN test (section 4).")
nb.md("## 2. Leakage between evidence and label, measured",
      "Each evidence family (dataset) alone, held out by orthogroup fold. A family that alone recovers "
      f"`C.STANDS_OUT` of what all evidence together recovers is flagged for a person to read: a "
      "restated label and one dominant biological signal look alike in the numbers. All Toxoplasma "
      "compartment labels come from hyperLOPIT (`compartment_source`), so none was predicted from "
      "sequence.")
nb.code("rows = []",
        "for org, t in (('Tg', 'compartment'), ('Tg', 'lopit_unified'), ('Pf', 'lopit_pf_location')):",
        "    fr = C.family_recovery(S.Context.shipped(org), t)",
        "    top = fr.iloc[0]",
        "    rows.append({'organism': org, 'label': t, 'all families': fr['all_families'].iloc[0],",
        "                 'strongest family': top['family'], 'its per-class recall': top['per_class'],",
        "                 'share of all': top['share_of_all'], 'flagged': ', '.join(fr[fr.stands_out].family) or 'none'})",
        "pd.DataFrame(rows).round(2)")
nb.md("## 3. Certainty that means what it says",
      "Support mapped to the held-out rate of being right (isotonic), with its calibration error "
      "measured out of fold. A recipe whose certainty is not calibrated makes no claims.")
nb.code("pd.DataFrame([{'strategy': k, 'n': (c := C.certainty_model(tg, k)).n, 'calibration error': c.calibration_error,",
        "               'calibrated': c.calibrated, 'certainty at support 0.3/0.6/0.9': tuple(np.round(c([.3, .6, .9]), 2))}",
        "              for k in ('feature_knn', 'supervised_classifier', 'random_forest', 'physical_partners')]).round(3)")
nb.md("Calibration held across tachyzoite expression quartiles and across distance from the labelled "
      "genes (measured 2026-10-04, ECE 0.03-0.05) -- but 42% of unlabelled genes on compartment lie "
      "farther from any labelled gene than 95% of labelled genes do. Claims for those carry no certainty "
      "('outside tested range').")
nb.md("## 4. Recipes: the generator independent evidence can test",
      "For each generator: the verifiers measured independent of it, the calibration error after "
      "verification, and on held-out genes how many claims reached 80% and how often they were right.")
nb.code("C.candidates('Tg', 'compartment').round(3)")
nb.md("## 5. The shipped claims",
      "Per label: the recipe, whether it is proven (calibrated end to end), and what it claims about "
      "genes with no label. Tested = an independent check agreed or disagreed.")
nb.code("rec = C.recipes()",
        "rec[['organism', 'target', 'generator', 'verifiers', 'proven', 'heldout_confident_precision', 'claims']].round(2)")
nb.code("cl = C.shipped()",
        "cl.groupby(['organism', 'target', 'status'], observed=True).size().unstack(fill_value=0)")
nb.md("Discoveries: tested, confidence at least 0.8 and at least twice the claimed class's prior.")
nb.code("d = pd.concat([C.discoveries(o) for o in ('Tg', 'Pf')])",
        "d.groupby(['organism', 'target'], observed=True).size()")
nb.code("ctx = S.Context.shipped('Tg')",
        "top = C.discoveries('Tg').head(15)",
        "top.assign(product=ctx.product(np.asarray([ctx.index[g.upper()] for g in top.gene_id.astype(str)])))[",
        "    ['gene_id', 'product', 'target', 'claim', 'certainty', 'confidence', 'prior', 'lift']].round(2)")
nb.write(os.path.join(ROOT, "notebooks", "claims_2026_10_04.ipynb"))
