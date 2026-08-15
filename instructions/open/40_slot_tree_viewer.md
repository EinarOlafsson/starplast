# 40 — A window that shows every slot, its place in the tree, and what fills it

**Status: open. Requested 2026-08-15. Execute after 39.**

## Why

The slot catalog is now 239 entries across three hierarchies, two organisms, and soon several
species and host tables. Every check run on it so far has been a Python one-liner written for the
occasion and thrown away — which is how the first audit of it reported a leak that did not exist,
because the check reimplemented pattern matching instead of calling `slots.declared_columns`.

**A sanity check nobody can see is a sanity check nobody performs.** This window is the audit made
permanent and visible: open it and the shape of the whole catalog is in front of you, including the
parts that are empty.

## What it is

`View → Slot tree`, opening a window of its own — not a tab, not a modal. Same rule as Preferences:
shown rather than exec'd, movable independently of the map, changes under it take effect while it
sits there. A modal here would mean checking the catalog against the map from memory.

### The tree

Clickable, one level per hierarchy level, leaves are slots:

    evidence ▾                                    [ evidence | biology | context ]
      molecular measurements ▾
        RNA ▾
          transcript abundance ▾
            Tg_transcription · tachyzoite      A   7,739  95.1%   3 columns
            Tg_transcription · bradyzoite      A   7,739  95.1%   11 columns
            Tg_transcription · merozoite       —       0      —   1 candidate
            Pf_transcription · gametocyte      —       0      —   no candidate

* The hierarchy selector switches between the three trees over the same slots. The same slot
  appears in all three at different addresses, which is the point of having three.
* A group row aggregates what is under it: how many slots, how many filled, how many genes covered.
* An organism filter, and a free-text filter that matches slot names, columns and candidate titles.

### Selecting a slot shows

* its `axis`, `context`, `unit`, `policy`, `role`, `target_family`
* its **columns in the loaded table**, from `slots.declared_columns` — never a re-implementation of
  the matching; this window must show what the code does, not what a copy of the code does
* its **candidates**: PMID, title, accession where there is one, and whether it is downloaded
* coverage and grade, computed live from the node table rather than read from the CSV
* its address in each of the three hierarchies

### The four things it must make impossible to miss

These are the audits currently done by hand. Put them where nobody has to run anything:

1. **Empty slots** — a distinct colour, counted at every group row. This is the map of what is not
   measured yet, and it is the most useful thing in the window.
2. **Columns claimed by two slots** — zero today; the window should show that count and turn red if
   it is ever not zero.
3. **Columns claimed by no slot** — likewise zero today.
4. **Candidates with no accession** — currently every Plasmodium candidate. A citation nobody can
   download is not a filled slot, and the window should not let that read as one.

## Implementation notes

* Read from `slots.all_slots()` and `slots.declared_columns()`. **Do not duplicate the catalog** —
  a second copy of the hierarchy in GUI code is a second thing to keep in step, and it will not be
  kept in step.
* `slots.relationship_tree(organism, hierarchy)` already builds the nested structure; the window is
  a `QTreeWidget` over its output plus a details pane.
* The window must open with no node table loaded, showing the catalog with coverage blank. It is
  also the tool for inspecting the *Plasmodium* arm, which has no table at all yet, and a viewer
  that requires data to show a slot cannot show the arm that most needs looking at.
* Tooltips on every column heading, per the project rule; the explanation of what a slot IS belongs
  here, since this is where somebody meets the concept for the first time.

## Tests

* The tree matches the catalog: every slot appears exactly once per hierarchy, and every group row's
  count equals the number of leaves beneath it.
* Coverage shown for a slot equals `slots.resolve` on the same node table — the window and the
  pipeline must not be able to disagree.
* The four alarm counts are computed from the same functions the audit uses, and a synthetic
  catalog with a deliberate double claim turns the count non-zero.
* Opens and populates with `nodes=None`.
* Headless construction, as everything else in this suite: no test may need a display.
* 100% coverage, no `pragma`.
