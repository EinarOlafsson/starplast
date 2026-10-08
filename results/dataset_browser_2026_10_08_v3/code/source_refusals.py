"""Recorded source refusals scoped to a question, not an entire publication."""

REFUSED_CANDIDATES = {
    # Eight GEO series proposed for the IFN-gamma macrophage slot, opened 2026-08-19. Not one is a
    # transcriptome of the parasite inside an activated macrophage: the term that matched is in the
    # summaries, not in the samples. The slot's verdict stands.
    ("Tg_transcription · in IFN-gamma macrophage", "GSE313582"):
        "GCN5b/PHD1 knockdown transcriptome in RH Ku80 tachyzoites, HFF host, no macrophage",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE313273"):
        "ChIP-seq of GCN5b, PHD1, MORC, HDAC3 and histone marks; not a transcriptome at all",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE313048"):
        "ATAC-seq of GCN5b knockdown tachyzoites in HFF monolayers",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE334992"):
        "RNA-seq of purified EXTRACELLULAR tachyzoites, wild type against one knockout, no host",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE329845"):
        "TgPRO knockout against wild type under actinomycin D; a redox-adaptation experiment",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE300509"):
        "ME49 grown on different host cell LINES; the host axis is cell line, not activation state",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE275112"):
        "AP2XII-8 knockdown and CUT&Tag; a cell-cycle transcription factor in HFF",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE266204"):
        "AP2XII-9 knockdown and CUT&Tag; likewise",
    # The six ribosome-profiling series proposed for the cell-cycle translation slot are the same
    # nine that the recorded sweep already opened one at a time. None is synchronised or sorted.
    ("Tg_translation · per cell-cycle phase", "GSE302108"):
        "5'UTR MPRA -- reporter constructs, not the endogenous translatome",
    ("Tg_translation · per cell-cycle phase", "GSE302107"):
        "ribosome profiling of unsynchronised tachyzoites; no cell-cycle axis",
    ("Tg_translation · per cell-cycle phase", "GSE243206"):
        "eIF4E1 depletion driving bradyzoite formation; a stage axis, not a cycle axis",
    ("Tg_translation · per cell-cycle phase", "GSE129869"):
        "host-context ribosome profiling, already shipped as its own dataset",
    ("Tg_translation · per cell-cycle phase", "GSE99395"):
        "intracellular against extracellular ribosome profiling, already shipped",
    ("Tg_translation · per cell-cycle phase", "GSE43722"):
        "the unfolded protein response; a stress axis",
    # Opened 2026-08-19 for `Pf_glycosylation`, and refused on three counts at once: the modification
    # is measured on RECOMBINANT protein, the substrate list is a domain prediction, and the biology
    # is the mosquito stages rather than the blood stage the slot names.
    ("Pf_glycosylation · asexual blood stage", "PXD033470"):
        "tryptophan C-mannosylation shown on recombinant SPATR and MTRAP, with the substrate set "
        "defined by TSR-domain content rather than measured in the parasite, and in the "
        "transmission stages rather than the asexual blood stage",
    # And the one proposed for protein turnover is the study that already fills a different slot.
    ("Tg_protein turnover", "PXD033642"):
        "the Ca2+-responsive thermal-shift proteome (PMID 35976251), which already fills "
        "`Tg_target engagement / thermal shift`; melting behaviour is not a degradation rate",
}


def records(catalog=None):
    """Return the existing scoped refusals with explicit organism and entity unit."""
    from . import slots
    catalog = slots.all_slots() if catalog is None else catalog
    questions = {slot.organism + '_' + slot.name: slot for slot in catalog}
    result = []
    for (question, accession), reason in sorted(REFUSED_CANDIDATES.items()):
        if question not in questions:
            continue
        slot = questions[question]
        result.append({'source_id': 'candidate:' + accession, 'organism': slot.organism,
            'unit': slot.unit, 'question': question, 'reason': reason})
    return result
