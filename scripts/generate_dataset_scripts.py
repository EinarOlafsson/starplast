#!/usr/bin/env python3
"""Write one script per dataset in the registry, into scripts/datasets/.

Generated rather than hand-written, because 28 hand-maintained files describing 28 registry entries
is 28 chances for a citation, a URL or a quirk to drift out of step with the registry that the README
and the methods are both generated from. Regenerate after changing `datasets.REGISTRY`:

    python scripts/generate_dataset_scripts.py

A test asserts that regenerating produces no diff, so a registry edit without a regeneration fails
the suite rather than leaving a stale script behind.

The per-dataset specifics that are NOT in the registry -- a tab-separated file with a .csv extension,
which sheet of a workbook holds the matrix, which module normalises the resulting columns -- live in
SPECIALS below, keyed by dataset. That keeps them in one reviewable place instead of scattered
through generated files where an edit would be overwritten on the next run.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from starplast import datasets  # noqa: E402

OUT_DIR = os.path.join(HERE, "datasets")

# Per-dataset handling that the registry does not carry.
#   sep    -- explicit separator, for files whose extension lies about their format
#   sheet  -- which sheet of a workbook holds the table
#   id_col -- the identifier column, where guessing by resolution rate picks the wrong one
#   by     -- the module that does the authoritative normalisation of this dataset's columns
SPECIALS: dict[str, dict] = {
    "xue_singlecell": {"sep": "\t", "by": "cellcycle.add_all()"},
    "stage_enriched": {"by": "cellcycle.add_all()"},
    # hyperLOPIT, measured and ortholog-transferred, all through the localisation module.
    "lopit_tgon": {"by": "localisation.lopit_labels()"},
    "lopit_pfal": {"by": "localisation.lopit_labels()"},
    "lopit_cpar": {"by": "localisation.lopit_labels()"},
    "gse108740": {"by": "expression.load_all()"},
    "gse206344": {"by": "expression.load_all()", "sheet": "1"},
    "proteome_pru": {"by": "screens.proteomics()"},
    "phosphosites": {"by": "screens.proteomics()"},
}
# Every CRISPR screen is normalised by the same module, so it is stated once rather than repeated.
for _k in ("crispr_invitro", "crispr_invivo_composite", "crispr_macrophage", "crispr_young2019",
           "gra17_synthlethal", "invivo_platform", "gra12", "hosttx_effectors"):
    SPECIALS.setdefault(_k, {})["by"] = "screens.crispr_screens()"

TEMPLATE = '''#!/usr/bin/env python3
"""{name}

{provides}

    level / kind : {level} / {kind}
    provides     : {columns}
    coverage     : {coverage}{citation}{pmid}{accession}{url}{path}{derived}
{note}
Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.{by}`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/{key}.py
"""
from _common import run

KEY = "{key}"

if __name__ == "__main__":
    run(KEY{args})
'''


def _wrap(text: str, width: int = 96, indent: str = "    ") -> str:
    """Wrap a note to the file's line length, preserving its own paragraph breaks."""
    import textwrap
    out = []
    for para in str(text).split("\n"):
        out.extend(textwrap.wrap(para, width=width - len(indent)) or [""])
    return "\n".join(indent + line for line in out)


def render(d) -> str:
    spec = SPECIALS.get(d.key, {})
    args = ""
    for name in ("sep", "sheet", "id_col"):
        if name in spec:
            args += f", {name}={spec[name]!r}"
    by = spec.get("by", "build_graph.load_nodes()")
    if by != "build_graph.load_nodes()":
        args += f", normalised_by={by!r}"
    note = f"\nQuirks that cost time once:\n{_wrap(d.note)}\n" if d.note else ""
    return TEMPLATE.format(
        key=d.key,
        name=d.name,
        provides=d.provides,
        level=d.level,
        kind=d.kind,
        columns=", ".join(d.columns) if d.columns else "(edges or build inputs only)",
        coverage=d.coverage,
        # Absent fields are omitted rather than rendered. Several entries genuinely have no citation
        # or no fixed path, and a line reading "citation : None" states something false about the
        # source instead of saying nothing about it.
        citation=f"\n    citation     : {d.citation}" if d.citation else "",
        pmid=f"\n    PMID         : {d.pmid}" if d.pmid else "",
        accession=f"\n    accession    : {d.accession}" if d.accession else "",
        url=f"\n    url          : {d.url}" if d.url else "",
        path=f"\n    local path   : {d.path}" if d.path else "",
        derived=(f"\n    derived from : {', '.join(d.derived_from)}" if d.derived_from else ""),
        note=note,
        by=by,
        args=args,
    )


def main(out_dir: str = OUT_DIR, log=print) -> list[str]:
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for d in datasets.REGISTRY:
        path = os.path.join(out_dir, f"{d.key}.py")
        text = render(d)
        # Only rewrite on a real change, so regenerating does not churn mtimes in git status.
        if not os.path.exists(path) or open(path).read() != text:
            open(path, "w").write(text)
        written.append(path)
    log(f"{len(written)} dataset scripts in {out_dir}")
    return written


if __name__ == "__main__":
    main()
