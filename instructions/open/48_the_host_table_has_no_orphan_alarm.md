# 48 — The orphan alarm never looks at the host table

**Status: open. Found 2026-08-19, while giving the host table its second and third columns.**

## What is wrong

`slot_tree.audit` counts four things, and one of them is the orphan: a shipped column no slot in
the catalogue claims. It walks `nodes` and only `nodes` --

```python
described = set(claimed) | {"gene_id"}
orphan = sum(1 for c in nodes.columns if c not in described)
```

-- so `host_proteins.parquet` has never been checked. It has shipped an unclaimed column the whole
time: `pv_enrichment_log2`, twelve host proteins from the vacuole-uptake study, answering no slot in
the catalogue. That is precisely the state the alarm exists to make impossible, and the alarm could
not see it because it was pointed at one table.

The parasite side found this class of bug once already: `protein stage share` was a measurement the
catalogue could not describe, and the orphan count is what caught it. The host side had no such
catch until now.

## What to do

1. **Extend the alarm.** Count host-table columns that no `host_gene` slot claims, adding them to
   the same `orphan` number the window already shows. Two things to get right:
   * the host slots of BOTH arms have to be consulted, not `all_slots(organism)` -- `rbc_*` is
     claimed only by a Plasmodium slot, and a per-organism walk would report it as an orphan in the
     Toxoplasma window;
   * `host_id` and `host_name` are the key and its label, so they are described by definition.
2. **Give `pv_enrichment_log2` a slot**, or decide it is bookkeeping and mark it so. It is a real
   question -- which host proteins are pulled to the parasitophorous vacuole -- and it is not one of
   the four families in `HOST_FAMILIES`, which are per-tissue and this is not.
3. **The denominator needs a fallback.** `host_row_space` keys off the tissue in the slot name
   (`"<family> · <tissue>"`), and a host slot that is not about a tissue has none. Fall back to the
   slot's own declared columns, which is what every other unit already does.
4. Regenerate the atlas: the slot count moves 271 -> 272, and anything asserting 271 moves with it.
