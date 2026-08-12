"""All four targets on the FULL proteome, same grid as the subsample runs.

The cell-cycle full-proteome run showed the subsample was flattering: mean F1 0.489 -> 0.396 and best
cluster purity 52% -> 31%. Smaller samples produce tighter, purer clusters, so every number computed on
a 3,000-gene subsample is optimistic -- including all four in the corrected table.

It also showed the cost objection was wrong. 90 full-proteome runs took 269 seconds, about 3s each,
which is what the subsample cost. There is no reason to quote subsample estimates in a paper when the
real thing is affordable.
"""
import json, sys, time
sys.path.insert(0, "/mnt/firecuda2/Claude/repo/starplast")
import pandas as pd
from starplast import paths, search
from starplast.tuning import EmbeddingStore

OUT = sys.argv[1]
nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
store = EmbeddingStore(f"{OUT}/embeddings")

TARGETS = [("cellcycle_phase", "measured"),
           ("compartment", "measured"),
           ("stage_enriched_derived", "POSITIVE control"),
           ("attention_depth", "NEGATIVE control")]

summary = {}
for target, role in TARGETS:
    print(f"\n{'='*78}\n{target}  --  {role}  (full proteome)\n{'='*78}", flush=True)
    t0 = time.time()
    res, per = search.search(
        nodes, target=target,
        na_policies=("median",), scalings=("rank",),
        n_neighbors_values=(15, 50), min_dist_values=(0.0, 0.25),
        min_cluster_sizes=(25, 60), sample_size=None, seed=42,
        store=store, save_above=0.0, log=print)
    res.to_csv(f"{OUT}/full_{target}.csv", index=False)
    if per is not None and len(per):
        per.to_csv(f"{OUT}/fulllab_{target}.csv", index=False)
    if len(res):
        b = res.sort_values("mean_f1", ascending=False).iloc[0]
        summary[target] = {"role": role, "mean_f1": round(float(b.mean_f1), 3),
                           "labels": int(b.n_labels_scored), "blocks": str(b.blocks),
                           "n_neighbors": int(b.n_neighbors), "min_dist": float(b.min_dist),
                           "min_cluster_size": int(b.min_cluster_size)}
        print(f"\n{target}: {b.mean_f1:.3f} over {int(b.n_labels_scored)} labels "
              f"in {time.time()-t0:.0f}s", flush=True)
    else:
        summary[target] = {"role": role, "mean_f1": None}
        print(f"\n{target}: no scorable runs", flush=True)

json.dump(summary, open(f"{OUT}/full_summary.json", "w"), indent=1)
SUB = {"cellcycle_phase": 0.489, "compartment": 0.228,
       "stage_enriched_derived": 0.709, "attention_depth": 0.339}
print("\n=== full proteome vs the 3,000-gene subsample ===")
print(f"{'target':24s} {'role':18s} {'subsample':>10s} {'full':>8s} {'change':>8s}")
for t, v in sorted(summary.items(), key=lambda kv: -(kv[1]["mean_f1"] or -1)):
    f = v["mean_f1"]
    ch = "" if f is None else f"{f - SUB[t]:+.3f}"
    print(f"{t:24s} {v['role']:18s} {SUB[t]:10.3f} {('-' if f is None else f'{f:.3f}'):>8s} {ch:>8s}")
